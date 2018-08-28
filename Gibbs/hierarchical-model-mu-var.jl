## Julia program for Hierarchical Modeling of means and variances
## author: weiya <szcfweiya@gmail.com>
## date: 27 August, 2018

using Distributions
using SpecialFunctions
using StatsBase
using DelimitedFiles

function higibbs(Y, T, mu0 = 50.0, gamma20 = 25.0, nu0 = 1.0, sigma20 = 100.0, eta0 = 1.0, tau20 = 100.0, a = 1.0, b = 100.0, alpha = 1.0, NUMAX = 40)
    m = size(unique(Y[:,1]), 1)
    # starting value
    ybar = ones(m)
    sv = ones(m)
    n = ones(m)
    for j = 1:m
        yj = Y[ [Y[i,1] == j for i = 1:end], 2]
        ybar[j] = mean(yj)
        sv[j] = var(yj)
        n[j] = size(yj, 1)
    end
    theta = ybar
    sigma2 = sv
#    sigma20 = 1 / mean(sigma2)
#    nu0 = 2 * mean(sigma2)^2 / var(sigma2)
    mu = mean(theta)
    tau2 = var(theta)

    THETA = ones(T, m)
    SIGMA2 = ones(T, m)
    # mu tau2 sigma20 nu0
    MTSN = ones(T, 4)

    for t = 1:T        
        # sample mu
        varmu = 1 / (m / tau2 + 1 / gamma20)
        meanmu = varmu * (m * mean(theta) / tau2 + mu0 / gamma20)
        rnorm = Normal(meanmu, sqrt(varmu))
        mu = rand(rnorm, 1)[1]

        # sample 1/tau2
        shapetau = (eta0 + m) / 2
        ratetau = ( eta0 * tau20 + sum((theta .- mu).^2) ) / 2
        rgamma = Gamma(shapetau, 1/ratetau)
        tau2 = 1 / rand(rgamma, 1)[1]

        # sample theta
        for j = 1:m
            vartheta = 1 / (n[j] / sigma2[j] + 1 / tau2)
            meantheta = vartheta * ( n[j] * mean(Y[ [Y[i,1] == j for i = 1:end], 2]) / sigma2[j] + mu / tau2) 
            rnorm = Normal(meantheta, sqrt(vartheta))
            theta[j] = rand(rnorm, 1)[1]
        end
        THETA[t, :] .= theta
        
        # sample sigma2
        for j = 1:m
            shapesig = (nu0 + n[j])/2
            yj = Y[ [Y[i,1] == j for i = 1:end], 2]
            ratesig = ( nu0*sigma20 + sum( (yj .- theta[j]).^2 ) )/2
            rgamma = Gamma(shapesig, 1/ratesig)
            sigma2[j] = 1 / rand(rgamma, 1)[1]            
        end
        SIGMA2[t, :] .= sigma2
        
        # sample sigma20
        shapesig = a + 0.5 * m * nu0
        ratesig = b + 0.5 * nu0 * sum(sigma2)
        rgamma = Gamma(shapesig, 1/ratesig)
        sigma20 = rand(rgamma, 1)[1]
        
        # sample nu0
        x = 1:NUMAX
        lpnu0 = ones(NUMAX)
        lpnu0 .= m * ( .5 * x .* log.(sigma20 * x / 2) .- lgamma.(x/2) ) .+ (x / 2 .- 1) * sum(log.(1 ./ sigma2)) .- x .* (alpha + .5 * sigma20 * sum(1 ./ sigma2))
        nu0 = sample(x, pweights(exp.(lpnu0 .- maximum(lpnu0))))
        
        # store results
        MTSN[t, :] .= [mu, tau2, sigma20, nu0]
    end
    return THETA, SIGMA2, MTSN
end

# run
Y = readdlm("math-score-Y.csv")
THETA, SIGMA2, MTSN = higibbs(Y, 100)

using PyPlot
plot(MTSN[:,3])
show()