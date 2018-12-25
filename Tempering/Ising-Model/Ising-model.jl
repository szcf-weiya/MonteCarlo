# ###################################################################################################
# Simulation studies for one-dimensional Ising model
# 
# author: weiya
# date: Dec 25, 2018
# 
# refer to Section 2.4.2 & Section 10.4 of
# Liu, J. S. (2008). Monte Carlo strategies in scientific computing (2. ed). New York, NY: Springer.
#
# ###################################################################################################
using Distributions
using StatsBase
using Plots
using LaTeXStrings

function exactMethod(M = 2000; beta = 10, d = 100)
    V = ones(d)
    V[1] = exp(beta) + exp(-beta)
    for t = 2:d
        V[t] = V[1] * V[t-1]
    end
    Z = 2*V[d]
    x = ones(Int, d)
    x[d] = 2*rand(Bernoulli(1/2))-1
    for t = (d-1):1
#        p1 = V[t] * exp(x[d] * 1)
#        p2 = V[t] * exp(x[d] * (-1))
        p1 = exp(x[d])
        p2 = exp(-x[d])
        x[t] = 2*rand(Bernoulli(p1/(p1+p2)))-1
    end
    return x
end

function pdfIsing(x, T)
    n = length(x)
    y = sum(x[1:(n-1)] .* x[2:n])
    return exp(y/T)
end

function logpdfIsing(x, T)
    n = length(x)
    y = sum(x[1:(n-1)] .* x[2:n])
    return y/T
end

function sampleX!(x::Array{Int}, T::Float64)
    d = length(x)
    for i = 1:d
        if i == 1
            p1 = exp( x[i+1] / T )
            p2 = exp( -x[i+1] / T )
        elseif i == d
            p1 = exp( x[i-1] / T )
            p2 = exp( -x[i-1] / T )
        else
            p1 = exp( (x[i-1] + x[i+1]) / T )
            p2 = exp( (-x[i-1] - x[i+1]) / T )
        end
        x[i] = 2*rand(Bernoulli( p1 / (p1+p2) )) - 1
    end
end

function sampleX!(x::Array{Int,2}, T::Array{Float64})
    I = length(T)
    d = size(x, 1)
    for j = 1:I
        for i = 1:d
            if i == 1
                p1 = exp( x[i+1,j] / T[j] )
                p2 = exp( -x[i+1,j] / T[j] )
            elseif i == d
                p1 = exp( x[i-1,j] / T[j] )
                p2 = exp( -x[i-1,j] / T[j] )
            else
                p1 = exp( (x[i-1,j] + x[i+1,j]) / T[j] )
                p2 = exp( (-x[i-1,j] - x[i+1,j]) / T[j] )
            end
            x[i,j] = 2*rand(Bernoulli( p1 / (p1+p2) )) - 1
        end
    end
end


function singleMCMC(M = 2000; T = 0.1, d = 100)
    # intial x
    x = 2 * rand(Bernoulli(1/2), d) .- 1
    for m = 1:M
        sampleX!(x, T)
    end
    return x
end
# result: always trap in one of the modes

function parallelTemper(M = 2000, T = [0.1, 0.2, 0.3, 0.4]; d = 100, alpha0 = T[1])
    I = length(T)
    # initial 
    x = 2*rand(Bernoulli(1/2), d, I) .- 1
    num1 = zeros(3)
    num2 = zeros(3)
    res = zeros(M, 4)
    for m = 1:M
        if rand() < alpha0 # parallel step
            #for i = 1:I
                # !!! Doesn't work !!!
                # y1 = x[:,i]
                # sampleX!(x[:,i], T[i])
                # y2 = x[:,i]
                # println(sum(abs.(y2.-y1)))
            #end
            sampleX!(x, T)
        else
#            idx = sample(1:I, 2, replace = false) # not neigbor
            idx1 = sample(1:(I-1))
            num1[idx1] += 1
            idx = [idx1, idx1+1]
#            rho = pdfIsing(x[:,idx[2]], T[idx[1]]) * pdfIsing(x[:,idx[1]], T[idx[2]]) / (pdfIsing(x[:,idx[1]], T[idx[1]]) * pdfIsing(x[:,idx[2]], T[idx[2]]))
            rho = logpdfIsing(x[:,idx[2]], T[idx[1]]) + logpdfIsing(x[:,idx[1]], T[idx[2]]) - (logpdfIsing(x[:,idx[1]], T[idx[1]]) + logpdfIsing(x[:,idx[2]], T[idx[2]]))
            if log(rand()) < rho
                # swappin step
                num2[idx1] += 1
                tmp = copy(x[:, idx[1]])
                x[:,idx[1]] .= x[:, idx[2]]
                x[:,idx[2]] .= tmp
            end
        end
        res[m, :] = sum(x, dims=1)
    end
    return res, num2 ./ num1
end

# run
for i in 1:10
    res, ratio = parallelTemper(200000)

    p1 = histogram(res[190000:end, 1], nbins = 20, title = L"T=0.1")
    p2 = histogram(res[190000:end, 2], nbins = 20, title = L"T=0.2")
    p3 = histogram(res[190000:end, 3], nbins = 20, title = L"T=0.3")
    p4 = histogram(res[190000:end, 4], nbins = 20, title = L"T=0.4")
    plot(p1, p2, p3, p4, layout = (2, 2),legend = false)
    savefig("IsingHist-rep$i.png")

    plot(autocor(res), label = [L"T=0.1", L"T=0.2", L"T=0.3", L"T=0.4"], title = "ACF")
    savefig("IsingACF-rep$i.png")
end 