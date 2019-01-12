using DelimitedFiles
using OffsetArrays
using SpecialFunctions
using StatsBase

data = readdlm("coal-mining-disaster.txt")
α = 1
β = 200/365.24
L = data[end] - data[1]
coal = data .- data[1]
λ = 3

# count number of disaster in [s1, s2)
function count_disasters(s1::Float64, s2::Float64)
    return(sum( (coal .>= s1) .& (coal .< s2) ))
end

# accept probability for H move
function alphaH(s::OffsetArray, h::OffsetArray, hjp::Float64, j::Int)
    m = count_disasters(s[j], s[j+1])
    log_ll_r = m * (log(hjp) - log(h[j])) - (hjp - h[j]) * (s[j+1] - s[j])
    log_prior_r = α * (log(hjp) - log(h[j])) - β * (hjp - h[j])
    return log_ll_r + log_prior_r
end

# accept probability for P move 
function alphaP(s::OffsetArray, h::OffsetArray, sjp::Float64, j::Int)
    m1 = count_disasters(s[j-1], s[j])
    m1p = count_disasters(s[j-1], sjp)
    m2 = count_disasters(s[j], s[j+1])
    m2p = count_disasters(sjp, s[j+1])
    log_ll_r = (m1p - m1) * log(h[j-1]) + (m2p - m2) * log(h[j]) - (sjp - s[j]) * (h[j-1] - h[j])
    log_prior_r = log(s[j+1] - sjp) + log(sjp - s[j-1]) - log(s[j+1] - s[j]) - log(s[j] - s[j-1])
    return log_ll_r + log_prior_r
end

# accept probability for birth move
function alphaB(s::OffsetArray, h::OffsetArray, sstar, h1p, h2p, j, k)
    m1p = count_disasters(s[j], sstar)
    m2p = count_disasters(sstar, s[j+1])
    mj = m1p + m2p
    log_ll_r = (m1p*log(h1p) + m2p*log(h2p) - (m1p + m2p)*log(h[j])) - (h1p*(sstar-s[j])+h2p*(s[j+1]-sstar)-h[j]*(s[j+1]-s[j]))
    log_prior_r = log(2(2k+3)λ) - 2log(L) + log(sstar-s[j]) + log(s[j+1]-sstar) - log(s[j+1]-s[j]) + α*log(β) - log(gamma(α)) + (α - 1) *(log(h1p) + log(h2p) - log(h[j])) - β*(h1p + h2p - h[j])
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
    log_prior_r = - (log(2(2k+3)λ) - 2log(L) + log(s[j+1]-s[j]) + log(s[j+2]-s[j+1]) - log(s[j+2]-s[j]) + α*log(β) - log(gamma(α)) + (α - 1) *(log(h[j]) + log(h[j+1]) - log(hjp)) - β*(h[j] + h[j+1] - hjp))
    log_proposal_r = -( log(death(k+1)) + log(L) - log(birth(k)) - log(k+1) )
    log_Jacobian = log(hjp) - 2 * log(h[j] + h[j+1])
    return log_ll_r + log_prior_r + log_proposal_r + log_Jacobian
end

function birth(k; c = 3.6/7)
    return c * min(1, λ/(k+1))
end

function death(k; c = 3.6/7)
    if k == 0
        return 0
    else
        return c * min(1, (k+1)/λ)
    end
end

function move(n::Int, h::OffsetArray, s::OffsetArray, k::Int)
    H = []
    S = []
    println(s)
    K = ones(n)
    K[1] = k
    push!(H, copy(h))
    push!(S, copy(s))
    for i = 2:n
        println(i)
        probP = 0
        if k > 0
            probP = 0.5 * (1 - birth(k) - death(k))
        end
        probH = 1 - birth(k) - death(k) - probP
        flag = rand()
        println("ok")
        if flag < probH
            println("ok1")
            j = sample(0:k)
            u = rand() - 0.5
            hjp = h[j] * exp(u)
            if log(rand()) < alphaH(s, h, hjp, j)
                h[j] = hjp
            end
        elseif flag < probH + probP
            println("ok2")
            j = sample(1:k)
            sjp = rand() * (s[j+1]-s[j-1]) + s[j-1]
            if log(rand()) < alphaP(s, h, sjp, j)
                s[j] = sjp
            end
        elseif flag < probH + probP + birth(k)
            println("ok3")
            sstar = rand() * L
            j = findfirst(s .>= sstar) - 1
            println(j)
            u = rand()
            h1p = h[j] * exp(-(s[j+1]-sstar)/(s[j+1]-s[j]) * log((1-u)/u))
            h2p = (1-u)/u*h1p
            println("ok31")
            if log(rand()) < alphaB(s, h, sstar, h1p, h2p, j, k)
                push!(s.parent, 0)
                s[(j+2):end] .= s[(j+1):(end-1)]
                s[j+1] = sstar
                push!(h.parent, 0)
                h[(j+2):end] .= h[(j+1):(end-1)]
                h[j+1] = h2p
                h[j] = h1p
                k = k + 1
            end
        else
            println("ok4")
            j = sample(1:k) - 1
            r = ( s[j+1] - s[j] ) / (s[j+2] - s[j] )
            hjp = h[j]^r * h[j+1]^(1-r)
            if log(rand()) < alphaD(s, h, hjp, j, k)
                k = k - 1
                h[j] = hjp
                # the first plus 1 for converting index
                # the seconde plus 1 for the removed one h_{j+1}
                deleteat!(h.parent, j+1+1)
                deleteat!(s.parent, j+1+1)
            end
        end
        K[i] = k
        push!(H, copy(h))
        push!(S, copy(s))
    end
    return K, H, S
end

K, H, S = move(44000, OffsetArray([1.7], -1), OffsetArray([0, L], -1), 0)

# results
using Plots
K = K[4001:end]
H = H[4001:end]
S = S[4001:end]
histogram(K, legend = false, xlab = "Number of change points")

# k = 1
S1 = vcat([parent(x)' for x in S[K.==1]]...)
S2 = vcat([parent(x)' for x in S[K.==2]]...)
S3 = vcat([parent(x)' for x in S[K.==3]]...)

using KernelDensity
y1 = kde(S1[:,2])
y2 = kde(vec(S2[:,2:3]))
y3 = kde(vec(S3[:,2:4]))
mm = extrema(union(S1[:,2], vec(S2[:,2:3]), vec(S3[:,2:4])))
p1 = plot(range(mm..., length=100), z->pdf(y1,z), title="k=1")
p2 = plot(range(mm..., length=100), z->pdf(y2,z), title="k=2")
p3 = plot(range(mm..., length=100), z->pdf(y3,z), title="k=3")
plot(p1, p2, p3, layout=(3,1), legend=false)

