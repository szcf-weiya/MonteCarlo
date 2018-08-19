## Julia program for Box-Muller algorithm
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-19

function boxmuller()
    x = ones(2)
    u = rand(2)
    logu = sqrt(-2*log(u[1]))
    x[1] = logu * cos(2*pi*u[2])
    x[2] = logu * sin(2*pi*u[2])
    return(x)
end

## example
boxmuller()