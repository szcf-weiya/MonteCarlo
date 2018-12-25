# ###################################################################################################
# Simulation studies of a bimodal distribution
# 
# author: weiya
# date: Dec 25, 2018
# 
# refer to Section 11.5.1 of
# Liu, J. S. (2008). Monte Carlo strategies in scientific computing (2. ed). New York, NY: Springer.
#
# ###################################################################################################
using LinearAlgebra
using Distributions
using StatsBase
using Plots
using LaTeXStrings


Σ1 = Matrix(1.0I, 2, 2)
Σ2 = [1 0.9; 0.9 1]
Σ3 = [1 -0.9; -0.9 1]
invΣ2 = [1 -0.9; -0.9 1] ./ 0.19
invΣ3 = [1 0.9; 0.9 1] ./ 0.19
μ1 = [0, 0]
μ2 = [-6, -6]
μ3 = [4, 4]
w = [0.34, 0.33, 0.33]

function targetpdf(x)
    model1 = MvNormal(μ1, Σ1)
    model2 = MvNormal(μ2, Σ2)
    model3 = MvNormal(μ3, Σ3)
    return pdf(model1, x) * w[1] + pdf(model2, x) * w[2] + pdf(model3, x) * w[3]
end

function mycontour()
    # refer to https://docs.juliaplots.org/latest/examples/gr/#contours
    x = -10:0.5:10
    y = -10:0.5:10
    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))
    Z = copy(X)
    for i in 1:size(X,1)
        for j in 1:size(X,2)
            Z[i,j] = targetpdf([X[i,j], Y[i,j]])
        end
    end
    contour(x, y, Z)
    savefig("contour.png")
end

function propose(x; a = 4)
    θ = rand() * 2π
    r = rand() * 4
    return [x[1] + r*cos(θ), x[2] + r*sin(θ)]
end

function iMH(M = 200000)
    # start points
    x = [rand() - 0.5, rand() - 0.5]
    X = ones(M, 2)
    num = 0
    for m = 1:M
        xs = propose(x)
        ρ = targetpdf(xs) / targetpdf(x)
        if rand() < ρ
            x = xs
            num += 1
        end
        X[m,:] = x
    end
    return X, num/M
end

function propose2(x; a = 4)
    r = rand() * 2a - a
    return x + r
end


function iMH2(M = 200000)
    # start points
    x1 = rand() - 0.5
    x2 = rand() - 0.5
    X = ones(M, 2)
    num = 0
    for m = 1:M
        x1s = propose2(x1)
        ρ = targetpdf([x1s, x2]) / targetpdf([x1, x2])
        if rand() < ρ
            x1 = x1s
            num += 1
        end
        x2s = propose2(x2)
        ρ = targetpdf([x1, x2s]) / targetpdf([x1, x2])
        if rand() < ρ
            x2 = x2s
        end
        X[m,:] = [x1, x2]
    end
    return X, num/M # why acceptance ratio is two times than iMH. Got it, because I miss the direction random selection.
end

# derivative of targetpdf
function dm(x; step = 0.05, nmaxstep = 100)
    term1 = w[1] / 2π * exp(-0.5*sum(x.^2))
    term2 = w[2] / (2π * 0.19) * exp( -0.5 * transpose(x.+6) * invΣ2 * (x.+6) )
    term3 = w[3] / (2π * 0.19) * exp( -0.5 * transpose(x.-4) * invΣ3 * (x.-4) )
    term = term1 .+ term2 * invΣ2 + term3 * invΣ3
    # negative of gradient direction (remove the minus symbol already)
    dx1 = term[1,1]*x[1] + term[1,2]*x[2]
    dx2 = term[2,1]*x[1] + term[2,2]*x[2]
    dx = [dx1, dx2]
    # find local maximum
    bst = targetpdf(x)
    i = 0
    while i < nmaxstep
        i += 1
        cur = targetpdf(x .- step * dx)
        if cur > bst
            bst = cur 
            x = x .- step * dx
        else
            break
        end
    end
    return x
end

function cgmc(M = 200000)
    # start points (each col is a stream)
    x = rand(2, 2) .- 0.5
    X = ones(M, 4)
    for m = 1:M
        # randomly choose xa
        xaIdx = sample([1,2])
        xa = x[:,xaIdx]
        # find gradient of pdf at xa (or conjugate gradient of log pdf at xa)
        y = dm(xa)
        xcIdx = 3 - xaIdx
        xc = x[:,xcIdx]
        e = y - xc
        e = e / norm(e)
        # sample r
        r, = MTM(y, e)
        # update xc
        x[:,xcIdx] = y .+ r*e
        X[m, :] .= vec(x)
    end
    return X
end

f(r, y, e) = abs(r) * targetpdf(y.+r*e)

function MTM(anchor, direction;M = 100, k = 5, σ = 10)
    # use the simplified version: OBMC algorithm
    # initial x
    x = 2*rand()-1
    num = 0
    for m = 1:M 
        y = ones(k)
        wy = ones(k)
        for i = 1:k
            y[i] = rand(Normal(x, σ))
            wy[i] = f(y[i], anchor, direction)
        end
        #wy = wy ./ sum(wy)
        yl = sample(y, pweights(wy))
        # draw reference points
        xprime = ones(k)
        xprime[k] = x
        wxprime = ones(k)
        wxprime[k] = f(x, anchor, direction)
        for i = 1:k-1
            xprime[i] = rand(Normal(yl, σ))
            wxprime[i] = f(xprime[i], anchor, direction)
        end
        # accept or not 
        ρ = sum(wy) / sum(wxprime)
        if rand() < ρ
            x = yl
            num += 1
        end
    end
    return x, num/M
end

# run
res, = iMH()
X = cgmc()
p1 = histogram(res[100000:end,1], xlab = L"x_1")
p2 = histogram(res[100000:end,2], xlab = L"x_2")
p3 = histogram(X[100000:end,1], xlab = L"x_1")
p4 = histogram(X[100000:end,2], xlab = L"x_2")

p5 = plot(autocor(res[100000:end,1]))
p6 = plot(autocor(X[100000:end,1]))
plot(p1, p3, p5, p6, layout = (2,2), legend = false)
