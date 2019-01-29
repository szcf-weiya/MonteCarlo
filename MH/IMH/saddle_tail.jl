using Distributions

p = 6
λ = 9

function K(t)
    return 2λ*t / (1-2t) - p / 2 * log(1-2t)
end

function D1K(t)
    return ( 4λ*t ) / (1-2t)^2 + ( 2λ + p ) / ( 1 - 2t )
end

function D2K(t)
    return 2(p*(1-2t) + 4λ) / (1-2t)^3
end

function logf(n, t)
    return n * (K(t) - t*D1K(t)) + 0.5log(D2K(t))
end

function logg(n, t)
    return -1n * D2K(0) * t^2/2
end

function mh_saddle(T::Int = 10000; n::Int = 1)
    Z = zeros(T)
    g = Normal(0, 1 / sqrt(n*D2K(0)))
    for t = 1:T-1
        z = rand(g)
        logr = logf(n, z) - logf(n, Z[t]) + logg(n, Z[t]) - logg(n, z)
        if log(rand()) < logr 
            Z[t+1] = z
        else
            Z[t+1] = Z[t]
        end
    end
    return Z
end

Z = mh_saddle(n = 1)

function tau(x)
    return ( -1p + 2x - sqrt(p^2 + 8λ * x) ) / (4x)
end

println(sum(Z .> tau(36.225)) / 10000)
println(sum(Z .> tau(40.542)) / 10000)
println(sum(Z .> tau(49.333)) / 10000)

# n = 10
Z = mh_saddle(n = 10)

println(sum(Z .> tau(113.6667)) / 10000)
println(sum(Z .> tau(102.063)) / 10000)
println(sum(Z .> tau(96.19335)) / 10000)


# n = 100
Z = mh_saddle(n = 100)

println(sum(Z .> tau(25.18054)) / 10000)
println(sum(Z .> tau(25.52361)) / 10000)
println(sum(Z .> tau(26.17395)) / 10000)
