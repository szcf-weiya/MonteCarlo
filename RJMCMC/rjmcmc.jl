using DelimitedFiles
using OffsetArrays
using SpecialFunctions

data = readdlm("coal-mining-disaster.txt")
α = 1
β = 200/365.24
L = data[end] - data[1]
coal = data .- data[1]
λ = 3

# count number of disaster in [s1, s2)
function count_disasters(s1::Float64, s2::Float64)
    return(sum(coal .>= s1 & coal .< s2))
end

# accept probability for H move
function alphaH(s::OffsetArray, j::Int, newh::Float64, h::Float64)
    m = count_disasters(s[j], s[j-1])
    log_ll_r = m * (log(newh) - log(h))
    log_prior_r = α * (log(newh) - log(h)) - β * (newh - h)
    return log_ll_r + log_prior_r
end

# accept probability for P move 
function alphaP(s::OffsetArray, j::Int, news::Float64, h::OffsetArray)
    m1 = count_disasters(s[j-1], s[j])
    m1p = count_disasters(s[j-1], news)
    m2 = count_disasters(s[j], s[j+1])
    m2p = count_disasters(news, s[j+1])
    log_ll_r = (m1p - m1) * log(h[j-1]) + (m2p - m2) * log(h[j]) - (news - s[j]) * (h[j-1] - h[j])
    log_prior_r = log(s[j+1] - news) + log(news - s[j-1]) - log(s[j+1] - s[j]) - log(s[j+1] - s[j])
    return log_ll_r + log_prior_r
end

# accept probability for birth move
function alphaB(s::OffsetArray, h::OffsetArray, sstar, h1p, h2p, j, k)
    m1p = count_disasters(s[j], sstar)
    m2p = count_disasters(sstar, s[j+1])
    mj = m1p + m2p
    log_ll_r = (m1p*log(h1p) + m2p*log(h2p) - (m1p + m2p)*log(h[j])) - (h1p*(sstar-s[j])+h2p*(s[j+1]-sstar)-h[j]*(s[j+1]-s[j]))
    log_prior_r = log(2(2k+3)λ) - 2log(L) + 2log(sstar-s[j]) + 2log(s[j+1]-sstar) - log(s[j+1]-s[j]) + α*log(β) - gamma(α) + (α - 1) *(log(h1p) + log(h2p) - log(h[j])) - β*(h1p + h2p - h[j])
    log_proposal_r = log(death(k+1)) + log(L) - log(birth(k)) - log(k+1)
    log_Jacobian = 2 * log(h1p + h2p) - log(h[j])
    return log_ll_r + log_prior_r + log_proposal_r + log_Jacobian
end

# accept probability for death move
function alphaD(s::OffsetArray, h::OffsetArray, hjp, j, k)
    m1 = count_disasters(s[j], s[j+1])
    m2 = count_disasters(s[j+1], s[j+2])
    mjp = m1 + m2
    log_ll_r = ( (m1+m2)*log(hjp) - m1*log(h[j]) - m2*log(h[j+1])) - (hjp*(s[j+2]-s[j]) - (h[j]*(s[j+1]-s[j])) - h[j+1]*(s[j+2]-s[j+1]))
    log_prior_r = - (log(2(2k+3)λ) - 2log(L) + 2log(s[j+1]-s[j]) + 2log(s[j+2]-s[j+1]) - log(s[j+2]-s[j]) + α*log(β) - gamma(α) + (α - 1) *(log(h[j]) + log(h[j+1]) - log(hjp)) - β*(h[j] + h[j+1] - hjp))
    log_proposal_r = -( log(death(k+1)) + log(L) - log(birth(k)) - log(k+1) )
    log_Jacobian = log(hjp) - 2 * log(h[j] + h[j+1])
    return log_ll_r + log_prior_r + log_proposal_r + log_Jacobian
end

function birth(k; c = 3.6/7)
    return c * min(1, λ/(k+1))
end

function death(k; c = 3.6/7)
    return c * min(1, (k+1)/λ)
end

