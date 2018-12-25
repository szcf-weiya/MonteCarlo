# #######################################
# Example 10.21 Probit Model of
# Robert, C. P., & Casella, G. (2010). Monte Carlo statistical methods (2. ed., softcover reprint of the hardcover 2. ed. 2004). New York, NY: Springer.
# #######################################

using Distributions
# generate data
function genD(n = 2000, sigma = sqrt(2.0), beta = 5.0)
    # sample r from Bernoulli(1/2)
    r = rand(Bernoulli(), n)
    d = ones(Int, n)
    for i = 1:n
        # sample z from N(r*beta, sigma2)
        zi = rand(Normal(r[i]*beta, sigma))
        if zi >= 0
            d[i] = 1
        else
            d[i] = 0
        end 
    end
    return r, d
end

# Gibbs sampler
function gibbs(R, D, M = 100)
    n = length(D)
    # initial
    beta = 25
    sigma = 5
    BETA = ones(M)
    SIGMA = ones(M)
    for m = 1:M
        # sample zi 
        z = ones(n)
        for i = 1:n
            z[i] = rand(TruncatedNormal(0, sigma / R[i], -Inf, beta)) * (2D[i] - 1)
        end
        # sample sigma
        rz = sum(R.*z)
        println(rz)
        invs2 = rand(Gamma(n/2+1.5, 1/( rz + 1/1.5 )) )
        sigma = 1 / sqrt(invs2)
        beta = rand(Normal(0, 1/( 2*(rz/sigma^2 + 1/200) )))
        SIGMA[m] = sigma 
        BETA[m] = beta
    end
    return BETA, SIGMA
end

# run 
R, D = genD()
gibbs(R, D)
