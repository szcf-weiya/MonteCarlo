using Distributions
a = 10
N = 1000
x = ones(N)
y = ones(N)
for t = 2:N
    u1 = rand() * exp(-abs(x[t-1]))
    u2 = rand() * exp(-abs(y[t-1]))
    u3 = rand() * exp(-a*abs(x[t-1]-y[t-1]))
    x[t] = rand(Uniform(log(u1), -log(u1)))
    while a*abs(x[t]-y[t-1]) > -log(u3)
        x[t] = rand(Uniform(log(u1), -log(u1)))
    end 
    y[t] = rand(Uniform(log(u2), -log(u2)))
    while a*abs(x[t]-y[t]) > -log(u3)
        y[t] = rand(Uniform(log(u2), -log(u2)))
    end
end
using Plots
plot(x, y, seriestype=:scatter, legend=false)