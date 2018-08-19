## Julia program for Atkinson's Poisson Simulation
## author: weiya <szcfweiya@gmail.com>
## date: 2018-08-19

function AtkinsonPois(lambda)
    # parameters
    beta = pi/sqrt(3*lambda)
    alpha = lambda*beta
    c = 0.767 - 3.36/lambda
    k = log(c) - lambda - log(beta)
    # step 1: propose new x
    u1 = rand()
    while true
        global x
        x = (alpha - log((1-u1)/u1))/beta
        x > -0.5 && break
    end
    while true
        # step 2: transform to N 
        N = floor(Int, x)
        u2 = rand()
        # step 3: accept or not
        lhs = alpha - beta*x + log(u2/(1+exp(alpha-beta*x))^2)
        rhs = k + N*log(lambda) - log(factorial(lambda))
        if lhs <= rhs
            return(N)
        end
    end
end

## example
N = 100;
res = ones(Int, N);
for i = 1:N
    res[i] = AtkinsonPois(10)
end
# ans: 8 9 13 10 12 .......


# alternative method
function SimplePois(lambda)
    s = 0
    k = 0
    while true
        u = rand()
        x = -log(u)/lambda
        s = s + x
        if s > 1
            return(k)
        end
        k = k + 1
    end
end

## example for simple poisson
res2 = ones(Int, N);
for i = 1:N 
    res2[i] = SimplePois(10)
end
# ans: 5 7 16 7 10 .......