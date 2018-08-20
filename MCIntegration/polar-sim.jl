## Julia program for Polar Simulation
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-20

using Statistics

function polarsim(x::Array)
    while true
        phi1 = rand()*2*pi
        phi2 = rand()*pi - pi/2
        u = rand()
        xdotxi = x[1]*cos(phi1) + x[2]*sin(phi1)*cos(phi2) + x[3]*sin(phi1)*sin(phi2)
        if u <= exp(xdotxi - sum(x.^2)/2)
            return(randn() + xdotxi)
        end
    end
end

# approximate E^pi
function Epi(m, x = [0.1, 1.2, -0.7])
    rhos = ones(m)
    for i = 1:m
        rhos[i] = 1/(2*polarsim(x)^2+3)
    end
    return(mean(rhos))
end

# example
Epi(10)