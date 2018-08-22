## Julia program for Adaptive Rejection Sampling 
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-20
## 
## refer to Inferentialist's R code: 
##   https://blog.inferentialist.com/2016/09/26/adaptive-sampling.html

module corears
# example 3.22 in Davison(2008)
r = 2; m = 10; mu = 0; sig2 = 1;
# range of y
yl = -2; yu = 2;
# function of h
function h(y)
    return(r * y - m * log(1 + exp(y)) - (y - mu)^2 / (2 * sig2))
end

function h(y::Array)
    return(r * y - m * log.(1 .+ exp.(y)) .- (y .- mu) .^2 ./ (2 * sig2))
end

# derivative of h
function dh(y)
    return(r - m * exp(y) / (1 + exp(y)) - (y - mu) / sig2)
end

function dh(y::Array)
    return(r .- m * exp.(y) ./ (1 .+ exp.(y)) .- (y .- mu)./sig2)
end

# intersection point
function zfix(yfixed::Array)
    yf0 = yfixed[1:end-1]
    yf1 = yfixed[2:end]
    zfixed = yf0 .+ (h(yf0) .- h(yf1) .+ (yf1 .- yf0) .* dh(yf1)) ./ (dh(yf1) .- dh(yf0))
    return(zfixed)
end

# evaluate log-density
function hplus(y, yfixed::Array)
    zfixed = zfix(yfixed)
    n = size(zfixed, 1)
    for i = 1:n
        if i == 1 && y < zfixed[i]
            return(h(yfixed[i]) + (y - yfixed[i]) * dh(yfixed[i]))
        elseif i < n && y >= zfixed[i] && y < zfixed[i+1]
            return(h(yfixed[i+1]) + (y - yfixed[i+1]) * dh(yfixed[i+1]))
        elseif i == n && y >= zfixed[n]
            return(h(yfixed[i+1]) + (y - yfixed[i+1]) * dh(yfixed[i+1]))
        end
    end
end

# calculate G_+(z_i)
function gplus_cdf(yfixed::Array, zfixed::Array)
    n = size(zfixed, 1)
    s = zeros(n+1)
    pr = zeros(n+1)
    for i = 1:(n+1)
        ## integral from -infty to zi
        if i == 1
        #    s[i] = exp(h(yfixed[i])) / dh(yfixed[i]) * (exp((zfixed[i]-yfixed[i]) * dh(yfixed[i])) - )
            s[i] = exp(h(yfixed[i])) / dh(yfixed[i]) * (exp((zfixed[i]-yfixed[i]) * dh(yfixed[i])) - exp((yl-yfixed[i]) * dh(yfixed[i])))
        elseif i == n+1
            s[i] = exp(h(yfixed[i])) / dh(yfixed[i]) * (exp((yu-yfixed[i]) * dh(yfixed[i])) - exp((zfixed[n]-yfixed[i]) * dh(yfixed[i])))
        else
            s[i] = exp(h(yfixed[i])) / dh(yfixed[i]) * (exp((zfixed[i]-yfixed[i]) * dh(yfixed[i])) - exp((zfixed[i]-yfixed[i]) * dh(yfixed[i])))
        end
    end
    pr = s / sum(s)
    return cumsum(pr), sum(s)
end

# sample from gplus density
function gplus_sample(yfixed::Array)
    zfixed = zfix(yfixed)
    gp = gplus_cdf(yfixed, zfixed)
    zpr = gp[1]
    norm_const = gp[2]
    n = size(zfixed, 1)
    u = rand()
    # Invert the gplus pdf
    for i = 1:n
        if i == 1 && u < zpr[i]
            ey = u * dh(yfixed[i]) * norm_const / exp(h(yfixed[i])) + exp((yl-yfixed[i])*dh(yfixed[i]))
            return(yfixed[i] + log(ey)/dh(yfixed[i]))    
        elseif i == n && u >= zpr[i]
            ey = (u - zpr[i]) * dh(yfixed[i+1]) * norm_const / exp(h(yfixed[i+1])) + exp((zfixed[i]-yfixed[i+1])*dh(yfixed[i+1]))
            return(yfixed[i+1] + log(ey)/dh(yfixed[i+1]))
        elseif i < n && u >= zpr[i] && u < zpr[i+1]
            ey = (u - zpr[i]) * dh(yfixed[i+1]) * norm_const / exp(h(yfixed[i+1])) + exp((zfixed[i]-yfixed[i+1])*dh(yfixed[i+1]))
            return(yfixed[i+1] + log(ey)/dh(yfixed[i+1]))
        end
    end
end


## adaptive rejection sampling
function ars(yfixed::Array)
    x = gplus_sample(yfixed)
    u = rand()
    if u <= exp(h(x)-hplus(x, yfixed))
        return(x)
    else
        return(ars(append!(yfixed, x)))
    end
end

## example
N = 100
res = ones(N);
for i = 1:N
    res[i] = ars([-1.8,-1.1,-0.5,-0.2])
end

export zfix, hplus, gplus_cdf, gplus_sample
end