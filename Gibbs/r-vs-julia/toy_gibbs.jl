function bigibbs(T::Int64, rho::Float64)
    x = ones(T+1)
    y = ones(T+1)
    for t = 1:T
        x[t+1] = randn() * sqrt(1-rho^2) + rho*y[t]
        y[t+1] = randn() * sqrt(1-rho^2) + rho*x[t+1]
    end
    return x, y
end

## example
x, y = bigibbs(Int64(2e6), 0.5)

# ## plot
# using PyPlot
# plt[:subplot](211)
# plt[:hist](x)
# ylabel(L"$x$")
# plt[:subplot](212)
# plt[:hist](y)
# ylabel(L"$y$")
# plt[:tight_layout]()
# plt[:savefig]("hist-gibbs-julia.png")