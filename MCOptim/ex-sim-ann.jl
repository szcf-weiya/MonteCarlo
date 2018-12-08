r = 0.5
function T(t)
    return 1/log(t)
end
# target function
function h(x)
    return (cos(50x) + sin(20x))^2
end

N = 2500
x = ones(N)
y = ones(N)
for t = 1:(N-1)
    # step 1
    at = max(x[t]-r, 0)
    bt = min(x[t]+r, 1)
    u = rand() * (bt - at) + at 
    # step 2
    rho = min(exp( (h(u) - h(x[t])) / T(t) ), 1)
    if rand() < rho
        x[t+1] = u
        y[t+1] = h(u)
    else
        x[t+1] = x[t]
        y[t+1] = y[t]
    end
end