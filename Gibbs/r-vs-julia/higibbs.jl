using Distributions
using SpecialFunctions
using StatsBase
using DelimitedFiles

function higibbs(Y, T, nu0 = 1, s20 = 100.0, eta0 = 1.0, t20 = 100.0, mu0 = 50.0, g20 = 25.0)
    m = size(unique(Y[:,1]), 1)
    # starting value
    ybar = ones(m)
    sv = ones(m)
    n = ones(Int, m)
    for j = 1:m
        yj = Y[ [Y[i,1] == j for i = 1:end], 2]
        ybar[j] = mean(yj)
        sv[j] = var(yj)
        n[j] = size(yj, 1)
    end
    theta = copy(ybar)
    sigma2 = mean(sv)
    mu = mean(theta)
    tau2 = var(theta)

    THETA = ones(T, m)
    # sigma2 mu tau2 
    SMT = ones(T, 3)

    for t = 1:T        
        # sample theta
        for j = 1:m
            vartheta = 1 / (n[j] / sigma2 + 1 / tau2)
            meantheta = vartheta * ( n[j] * ybar[j] / sigma2 + mu / tau2) 
            rnorm = Normal(meantheta, sqrt(vartheta))
            theta[j] = rand(rnorm, 1)[1]
        end
        
        # sample new sigma2
        nun = nu0 + sum(n)
        ss = nu0 * s20
        for j = 1:m
            ss += sum((Y[ [Y[i, 1] == j for i = 1:end], 2] .- theta[j]).^2)
        end  
        rgamma = Gamma(nun/2, 2/ss)
        sigma2 = 1 / rand(rgamma, 1)[1]

        # sample mu
        varmu = 1 / (m / tau2 + 1 / g20)
        meanmu = varmu * (m * mean(theta) / tau2 + mu0 / g20)
        rnorm = Normal(meanmu, sqrt(varmu))
        mu = rand(rnorm, 1)[1]

        # sample 1/tau2
        shapetau = (eta0 + m) / 2
        ratetau = ( eta0 * t20 + sum((theta .- mu).^2) ) / 2
        rgamma = Gamma(shapetau, 1/ratetau)
        tau2 = 1 / rand(rgamma, 1)[1]
                
        # store results
        THETA[t, :] .= theta
        SMT[t, :] .= [sigma2, mu, tau2]
    end
    return THETA, SMT
end

# run
Y = readdlm("math-score-Y.csv")
THETA, SMT = higibbs(Y, 50);
#using PyPlot
#plt[:hist](SMT[:,2])
#show()